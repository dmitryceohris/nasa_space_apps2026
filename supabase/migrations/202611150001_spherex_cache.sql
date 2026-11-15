create schema if not exists spherex;
create table if not exists spherex.results (
    key text primary key,
    version integer not null,
    payload jsonb not null,
    expires_at timestamptz not null
);
create table if not exists spherex.images (
    result_key text references spherex.results(key) on delete cascade,
    filename text not null,
    png bytea not null,
    primary key (result_key, filename)
);
